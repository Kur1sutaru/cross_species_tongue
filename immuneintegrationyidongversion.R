setwd("/projects/r_workspace/cristal_data/cross_species/subtypes_cross_species/immune")
# =============================================================================
# Immune Compartment Re-integration & Subtype Annotation
# Input: existing `immune` Seurat v5 object (4 species)
# Strategy: Human + Marmoset → Harmony | Mouse + Rat → Harmony → annotate
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(scales)
})

set.seed(42)

# ─── Settings ─────────────────────────────────────────────────────────────────
outdir <- "immune_reintegration"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

RESOLUTION_PRIMATE <- 0.5
RESOLUTION_RODENT  <- 0.4
N_PCS              <- 30
N_HARMONY_DIMS     <- 25

species_colors <- c(
  Human    = "#1D9E75",
  Marmoset = "#377EB8",
  Mouse    = "#E41A1C",
  Rat      = "#984EA3"
)

theme_pub <- function(base_size = 11) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line         = element_line(linewidth = 0.4),
      axis.ticks        = element_line(linewidth = 0.3),
      axis.text         = element_text(color = "black"),
      legend.background = element_blank(),
      strip.background  = element_blank(),
      strip.text        = element_text(face = "bold"),
      plot.title        = element_text(face = "bold", hjust = 0.5)
    )
}

# =============================================================================
# STEP 0 — Confirm object
# =============================================================================
cat("── Step 0: Confirming immune object ──\n")
cat("Total cells:", ncol(immune), "\n")
cat("Species breakdown:\n")
print(table(immune$species))
cat("Metadata columns:\n")
print(colnames(immune@meta.data))

# =============================================================================
# STEP 1 — Split into primate and rodent
# =============================================================================
cat("\n── Step 1: Splitting into primate and rodent groups ──\n")

immune_primate <- subset(immune, subset = species %in% c("Human", "Marmoset"))
immune_rodent  <- subset(immune, subset = species %in% c("Mouse",  "Rat"))

cat("Primate — Human:", sum(immune_primate$species == "Human"),
    "| Marmoset:", sum(immune_primate$species == "Marmoset"),
    "| Total:", ncol(immune_primate), "\n")
cat("Rodent  — Mouse:", sum(immune_rodent$species == "Mouse"),
    "| Rat:", sum(immune_rodent$species == "Rat"),
    "| Total:", ncol(immune_rodent), "\n")

# =============================================================================
# STEP 2 — Pre-processing helper (NormalizeData → PCA)
# =============================================================================
preprocess_seurat <- function(seu, group_name, n_pcs = N_PCS) {
  cat(sprintf("\n  [%s] Normalizing...\n", group_name))
  seu <- NormalizeData(seu, verbose = FALSE)
  
  cat(sprintf("  [%s] Finding variable features (split by species)...\n", group_name))
  seu[["RNA"]] <- split(seu[["RNA"]], f = seu$species)
  seu <- FindVariableFeatures(seu, selection.method = "vst",
                              nfeatures = 3000, verbose = FALSE)
  
  cat(sprintf("  [%s] Scaling...\n", group_name))
  seu <- ScaleData(seu, verbose = FALSE)
  
  cat(sprintf("  [%s] PCA...\n", group_name))
  seu <- RunPCA(seu, npcs = n_pcs, verbose = FALSE)
  
  return(seu)
}

# =============================================================================
# STEP 3 — Harmony integration helper
# =============================================================================
run_harmony <- function(seu, group_name, resolution,
                        n_harmony_dims = N_HARMONY_DIMS) {
  cat(sprintf("\n  [%s] Running Harmony...\n", group_name))
  seu <- IntegrateLayers(
    object         = seu,
    method         = HarmonyIntegration,
    orig.reduction = "pca",
    new.reduction  = "harmony",
    group.by.vars  = "species",
    verbose        = TRUE
  )
  
  cat(sprintf("  [%s] Joining layers...\n", group_name))
  seu[["RNA"]] <- JoinLayers(seu[["RNA"]])
  
  cat(sprintf("  [%s] FindNeighbors...\n", group_name))
  seu <- FindNeighbors(seu, reduction = "harmony",
                       dims = 1:n_harmony_dims, verbose = FALSE)
  
  cat(sprintf("  [%s] FindClusters (res = %s)...\n", group_name, resolution))
  seu <- FindClusters(seu, resolution = resolution, verbose = FALSE)
  
  cat(sprintf("  [%s] UMAP...\n", group_name))
  seu <- RunUMAP(seu, reduction = "harmony",
                 dims = 1:n_harmony_dims, verbose = FALSE)
  
  cat(sprintf("  [%s] → %d clusters\n",
              group_name, length(levels(seu$seurat_clusters))))
  print(table(seu$seurat_clusters))
  return(seu)
}

# ── Run pipeline ──────────────────────────────────────────────────────────────
immune_primate <- preprocess_seurat(immune_primate, "Primate")
immune_primate <- run_harmony(immune_primate, "Primate", RESOLUTION_PRIMATE)
saveRDS(immune_primate, file.path(outdir, "immune_primate.rds"))

immune_rodent  <- preprocess_seurat(immune_rodent,  "Rodent")
immune_rodent  <- run_harmony(immune_rodent,  "Rodent",  RESOLUTION_RODENT)
saveRDS(immune_rodent, file.path(outdir, "immune_rodent.rds"))

cat("\n  Both objects integrated and saved.\n")

# =============================================================================
# STEP 4 — Canonical marker gene list (PanglaoDB / Tabula Muris / Tabula Sapiens)
# =============================================================================
marker_genes <- list(
  "Pan T"                  = c("CD3D", "CD3E", "CD3G"),
  "CD4+ T"                 = c("CD4", "IL7R", "TCF7"),
  "CD8+ T"                 = c("CD8A", "CD8B"),
  "Treg"                   = c("FOXP3", "IL2RA", "CTLA4"),
  "Naive T"                = c("SELL", "CCR7", "TCF7"),
  "Effector/Cytotoxic T"   = c("GZMB", "GZMK", "PRF1"),
  "Exhausted T"            = c("PDCD1", "LAG3", "HAVCR2"),
  "gdT"                    = c("TRGC1", "TRDV1"),
  "NKT"                    = c("KLRB1"),
  "Th17"                   = c("RORC"),
  "Tfh"                    = c("CXCR5", "BCL6"),
  "Pan B"                  = c("CD19", "MS4A1", "CD79A", "CD79B"),
  "Naive B"                = c("IGHD", "IGHM", "TCL1A"),
  "Memory B"               = c("CD27", "FCER2"),
  "Plasma cell"            = c("MZB1", "SDC1", "PRDM1", "XBP1", "IGHG1"),
  "NK"                     = c("NCAM1", "NKG7", "GNLY", "KLRD1", "KLRF1",
                               "FCGR3A", "TYROBP"),
  "Classical Mono"         = c("CD14", "LYZ", "S100A8", "S100A9", "VCAN"),
  "Non-classical Mono"     = c("FCGR3A", "CX3CR1"),
  "Pan Macro"              = c("CD68", "MRC1"),
  "Tissue-resident Macro"  = c("MARCO", "C1QA", "C1QB", "FOLR2"),
  "LAM Macro"              = c("TREM2", "APOE"),
  "M1 Macro"               = c("TNF", "IL1B"),
  "cDC1"                   = c("CLEC9A", "XCR1"),
  "cDC2"                   = c("CD1C", "CLEC10A", "FCER1A"),
  "pDC"                    = c("LILRA4", "CLEC4C", "IL3RA", "SIGLEC6"),
  "Migratory DC"           = c("CCR7", "CCL17", "CD83"),
  "Pan DC"                 = c("ITGAX", "HLA-DRA"),
  "Mast"                   = c("TPSAB1", "TPSB2", "KIT", "CPA3",
                               "HDC", "GATA2", "MS4A2", "FCER1G"),
  "Neutrophil"             = c("CXCR1", "CXCR2", "CSF3R", "FCGR3B",
                               "MPO", "ELANE"),
  "Pan leukocyte"          = c("PTPRC")
)

# =============================================================================
# STEP 5 — Module scores per lineage
# =============================================================================
cat("\n── Step 5: Module scores ──\n")

add_module_scores <- function(seu, markers, obj_name) {
  for (nm in names(markers)) {
    genes_present <- intersect(markers[[nm]], rownames(seu))
    if (length(genes_present) == 0) next
    safe_nm <- gsub("[^A-Za-z0-9_]", "_", nm)
    seu <- AddModuleScore(seu,
                          features = list(genes_present),
                          name     = paste0("score_", safe_nm),
                          ctrl     = 50, seed = 42)
    colnames(seu@meta.data)[ncol(seu@meta.data)] <- paste0("score_", safe_nm)
  }
  return(seu)
}

immune_primate <- add_module_scores(immune_primate, marker_genes, "Primate")
immune_rodent  <- add_module_scores(immune_rodent,  marker_genes, "Rodent")

# =============================================================================
# STEP 6 — FindAllMarkers (save CSVs for inspection)
# =============================================================================
cat("\n── Step 6: FindAllMarkers ──\n")

Idents(immune_primate) <- "seurat_clusters"
markers_primate <- FindAllMarkers(immune_primate, only.pos = TRUE,
                                  min.pct = 0.25, logfc.threshold = 0.3,
                                  test.use = "wilcox", verbose = FALSE)
write.csv(markers_primate,
          file.path(outdir, "markers_primate_clusters.csv"), row.names = FALSE)

Idents(immune_rodent) <- "seurat_clusters"
markers_rodent  <- FindAllMarkers(immune_rodent,  only.pos = TRUE,
                                  min.pct = 0.25, logfc.threshold = 0.3,
                                  test.use = "wilcox", verbose = FALSE)
write.csv(markers_rodent,
          file.path(outdir, "markers_rodent_clusters.csv"),  row.names = FALSE)

cat("  Marker CSVs saved to", outdir, "\n")

# =============================================================================
# STEP 7 — Exploratory UMAPs + DotPlot + FeaturePlot
# =============================================================================
cat("\n── Step 7: Exploratory plots ──\n")

dotplot_genes <- c(
  "CD3E", "CD4", "CD8A", "FOXP3", "GZMB", "PDCD1", "TRGC1", "KLRB1",
  "RORC", "BCL6", "CD19", "MS4A1", "IGHD", "CD27", "MZB1", "SDC1",
  "NCAM1", "NKG7", "GNLY", "CD14", "LYZ", "S100A8", "CD68", "MRC1",
  "C1QA", "TREM2", "IL1B", "ITGAX", "CLEC9A", "CD1C", "LILRA4",
  "HLA-DRA", "TPSAB1", "KIT", "CPA3", "MPO", "CSF3R", "FCGR3B", "PTPRC"
)

key_features <- c("CD3E", "CD19", "NKG7", "CD14",
                  "CD68", "ITGAX", "TPSAB1", "MPO", "PTPRC")

make_plots <- function(seu, group_name, spc_colors, dot_genes, feat_genes) {
  
  # UMAP by cluster + species
  p_clust <- DimPlot(seu, group.by = "seurat_clusters",
                     label = TRUE, label.size = 3.5, repel = TRUE) +
    ggtitle(paste(group_name, "— clusters")) + theme_pub() +
    theme(legend.position = "none")
  p_spc   <- DimPlot(seu, group.by = "species", cols = spc_colors) +
    ggtitle("Species") + theme_pub()
  
  pdf(file.path(outdir, paste0("01_umap_", tolower(group_name), "_exploratory.pdf")),
      width = 14, height = 6)
  print(p_clust | p_spc)
  dev.off()
  
  # DotPlot
  genes_ok <- intersect(dot_genes, rownames(seu))
  Idents(seu) <- "seurat_clusters"
  p_dot <- DotPlot(seu, features = genes_ok, dot.scale = 5) +
    RotatedAxis() +
    scale_color_gradient2(low = "#313695", mid = "white", high = "#a50026",
                          midpoint = 0) +
    ggtitle(paste(group_name, "— canonical markers")) +
    theme_pub(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          legend.key.size = unit(0.35, "cm"))
  
  pdf(file.path(outdir, paste0("02_dotplot_", tolower(group_name), ".pdf")),
      width = 18, height = 7)
  print(p_dot)
  dev.off()
  
  # FeaturePlot
  feats_ok <- intersect(feat_genes, rownames(seu))
  p_feat <- FeaturePlot(seu, features = feats_ok,
                        order = TRUE, ncol = 3,
                        cols = c("lightgrey", "#a50026"), pt.size = 0.3) &
    theme_pub(base_size = 9)
  
  pdf(file.path(outdir, paste0("03_featureplot_", tolower(group_name), ".pdf")),
      width = 14, height = 14)
  print(p_feat)
  dev.off()
  
  cat(sprintf("  [%s] Exploratory plots saved.\n", group_name))
}

make_plots(immune_primate, "Primate", species_colors, dotplot_genes, key_features)
make_plots(immune_rodent,  "Rodent",  species_colors, dotplot_genes, key_features)

# =============================================================================
# STEP 8 — ANNOTATION
# ─────────────────────────────────────────────────────────────────────────────
# INSTRUCTIONS:
#   1. Open the CSVs from Step 6 and the PDFs from Step 7
#   2. For each cluster, identify the top marker genes and cross-reference
#      with the marker_genes list above
#   3. Fill in the named vectors below (one entry per cluster number)
#   4. Re-run from here down
# =============================================================================
cat("\n── Step 8: Annotation (fill vectors below after marker inspection) ──\n")

# Check how many clusters you have in each object:
cat("Primate clusters:", levels(immune_primate$seurat_clusters), "\n")
cat("Rodent clusters: ", levels(immune_rodent$seurat_clusters),  "\n")

# ─── Fill these in ────────────────────────────────────────────────────────────
primate_annotations <- c(
  # "0" = "CD4+ T cell",
  # "1" = "CD8+ T cell",
  # "2" = "B cell",
  # ... add one line per cluster
)

rodent_annotations <- c(
  # "0" = "CD4+ T cell",
  # "1" = "Macrophage",
  # ... add one line per cluster
)
# ─────────────────────────────────────────────────────────────────────────────

# Apply annotations (run after filling vectors above)
if (length(primate_annotations) > 0) {
  immune_primate$immune_subtype <- unname(
    primate_annotations[as.character(immune_primate$seurat_clusters)]
  )
  cat("  Primate annotations applied.\n")
} else {
  immune_primate$immune_subtype <- paste0("Cluster_",
                                          immune_primate$seurat_clusters)
  cat("  Primate: placeholder subtype labels set (fill annotations above).\n")
}

if (length(rodent_annotations) > 0) {
  immune_rodent$immune_subtype <- unname(
    rodent_annotations[as.character(immune_rodent$seurat_clusters)]
  )
  cat("  Rodent annotations applied.\n")
} else {
  immune_rodent$immune_subtype <- paste0("Cluster_",
                                         immune_rodent$seurat_clusters)
  cat("  Rodent: placeholder subtype labels set (fill annotations above).\n")
}

# ─── Broad lineage ────────────────────────────────────────────────────────────
lineage_map <- c(
  "CD4+ T cell"                = "T/NK lymphoid",
  "CD8+ T cell"                = "T/NK lymphoid",
  "Treg"                       = "T/NK lymphoid",
  "Naive T cell"               = "T/NK lymphoid",
  "Effector T cell"            = "T/NK lymphoid",
  "Exhausted T cell"           = "T/NK lymphoid",
  "gdT cell"                   = "T/NK lymphoid",
  "NKT cell"                   = "T/NK lymphoid",
  "Th17"                       = "T/NK lymphoid",
  "Tfh"                        = "T/NK lymphoid",
  "NK cell"                    = "T/NK lymphoid",
  "B cell"                     = "B/Plasma",
  "Naive B cell"               = "B/Plasma",
  "Memory B cell"              = "B/Plasma",
  "Plasma cell"                = "B/Plasma",
  "Classical monocyte"         = "Myeloid",
  "Non-classical monocyte"     = "Myeloid",
  "Macrophage"                 = "Myeloid",
  "Tissue-resident macrophage" = "Myeloid",
  "LAM macrophage"             = "Myeloid",
  "M1 macrophage"              = "Myeloid",
  "cDC1"                       = "Myeloid",
  "cDC2"                       = "Myeloid",
  "pDC"                        = "Myeloid",
  "Migratory DC"               = "Myeloid",
  "Neutrophil"                 = "Myeloid",
  "Mast cell"                  = "Mast cell"
)

assign_lineage <- function(subtypes, lmap) {
  lin <- lmap[subtypes]
  lin[is.na(lin)] <- "Unannotated"
  unname(lin)
}

immune_primate$immune_lineage <- assign_lineage(immune_primate$immune_subtype, lineage_map)
immune_rodent$immune_lineage  <- assign_lineage(immune_rodent$immune_subtype,  lineage_map)

# =============================================================================
# STEP 9 — Annotated UMAPs
# =============================================================================
cat("\n── Step 9: Annotated UMAPs ──\n")

subtype_palette <- c(
  "CD4+ T cell"                = "#4DAF4A",
  "CD8+ T cell"                = "#377EB8",
  "Treg"                       = "#A65628",
  "Naive T cell"               = "#80B1D3",
  "Effector T cell"            = "#1D9E75",
  "Exhausted T cell"           = "#BEBADA",
  "gdT cell"                   = "#FDB462",
  "NKT cell"                   = "#B3DE69",
  "Th17"                       = "#FCCDE5",
  "Tfh"                        = "#D9D9D9",
  "NK cell"                    = "#BC80BD",
  "B cell"                     = "#E41A1C",
  "Naive B cell"               = "#FB8072",
  "Memory B cell"              = "#FFED6F",
  "Plasma cell"                = "#FF7F00",
  "Classical monocyte"         = "#984EA3",
  "Non-classical monocyte"     = "#CAB2D6",
  "Macrophage"                 = "#6A3D9A",
  "Tissue-resident macrophage" = "#33A02C",
  "LAM macrophage"             = "#B2DF8A",
  "M1 macrophage"              = "#A6CEE3",
  "cDC1"                       = "#1F78B4",
  "cDC2"                       = "#FDBF6F",
  "pDC"                        = "#FF6666",
  "Migratory DC"               = "#E31A1C",
  "Neutrophil"                 = "#B15928",
  "Mast cell"                  = "#00BFC4",
  "Unannotated"                = "#CCCCCC"
)

plot_annotated <- function(seu, group_name, palette, spc_colors) {
  p1 <- DimPlot(seu, group.by = "immune_subtype", cols = palette,
                label = TRUE, label.size = 3, repel = TRUE) +
    ggtitle(paste(group_name, "— immune subtypes")) + theme_pub() +
    theme(legend.text = element_text(size = 8),
          legend.key.size = unit(0.4, "cm"))
  p2 <- DimPlot(seu, group.by = "species", cols = spc_colors) +
    ggtitle("Species") + theme_pub()
  p3 <- DimPlot(seu, group.by = "immune_lineage",
                label = TRUE, label.size = 3.5, repel = TRUE) +
    ggtitle("Lineage") + theme_pub()
  p1 | p2 | p3
}

pdf(file.path(outdir, "04_umap_primate_annotated.pdf"), width = 20, height = 6)
print(plot_annotated(immune_primate, "Primate (Human + Marmoset)",
                     subtype_palette, species_colors))
dev.off()

pdf(file.path(outdir, "04_umap_rodent_annotated.pdf"), width = 20, height = 6)
print(plot_annotated(immune_rodent, "Rodent (Mouse + Rat)",
                     subtype_palette, species_colors))
dev.off()

cat("  Annotated UMAPs saved.\n")

# =============================================================================
# STEP 10 — Propagate annotations back to full immune object
# =============================================================================
cat("\n── Step 10: Propagating annotations to full immune object ──\n")

all_meta <- rbind(
  data.frame(cell           = colnames(immune_primate),
             immune_subtype = immune_primate$immune_subtype,
             immune_lineage = immune_primate$immune_lineage,
             stringsAsFactors = FALSE),
  data.frame(cell           = colnames(immune_rodent),
             immune_subtype = immune_rodent$immune_subtype,
             immune_lineage = immune_rodent$immune_lineage,
             stringsAsFactors = FALSE)
)

idx <- match(colnames(immune), all_meta$cell)
immune$immune_subtype <- all_meta$immune_subtype[idx]
immune$immune_lineage <- all_meta$immune_lineage[idx]

cat("Annotation check:\n")
print(table(immune$immune_subtype, useNA = "ifany"))

# =============================================================================
# STEP 11 — Save all objects
# =============================================================================
cat("\n── Step 11: Saving objects ──\n")
saveRDS(immune_primate, file.path(outdir, "immune_primate_annotated.rds"))
saveRDS(immune_rodent,  file.path(outdir, "immune_rodent_annotated.rds"))
saveRDS(immune,         file.path(outdir, "immune_all_annotated.rds"))
cat("  All objects saved to:", outdir, "\n")

cat("\n── Session info ──\n")
sessionInfo()



# =============================================================================
# Integrate immune_primate + immune_rodent into a single object
# =============================================================================

# ── Step 1: Merge ─────────────────────────────────────────────────────────────
cat("Merging primate and rodent immune objects...\n")

immune <- merge(
  x          = immune_primate,
  y          = immune_rodent,
  merge.data = TRUE
)

cat("Total cells after merge:", ncol(immune), "\n")
print(table(immune$species))

# ── Step 2: Re-process on merged object ───────────────────────────────────────
cat("Normalizing...\n")
immune <- NormalizeData(immune, verbose = FALSE)

cat("Finding variable features (split by species)...\n")
immune[["RNA"]] <- split(immune[["RNA"]],
                                    f = immune$species)
immune <- FindVariableFeatures(immune,
                                          selection.method = "vst",
                                          nfeatures = 3000,
                                          verbose = FALSE)

cat("Scaling...\n")
immune <- ScaleData(immune, verbose = FALSE)

cat("PCA...\n")
immune <- RunPCA(immune, npcs = 30, verbose = FALSE)

# ── Step 3: Harmony — correct for species (all 4) ─────────────────────────────
# group.by.vars = "species" handles all four species as batch variables
cat("Running Harmony...\n")
immune <- IntegrateLayers(
  object         = immune,
  method         = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction  = "harmony",
  group.by.vars  = "species",
  verbose        = TRUE
)

cat("Joining layers...\n")
immune[["RNA"]] <- JoinLayers(immune[["RNA"]])

# ── Step 4: Cluster + UMAP ────────────────────────────────────────────────────
cat("FindNeighbors...\n")
immune <- FindNeighbors(immune,
                                   reduction = "harmony",
                                   dims = 1:25,
                                   verbose = FALSE)

cat("FindClusters...\n")
immune <- FindClusters(immune,
                                  resolution = 0.5,
                                  verbose = FALSE)

cat("UMAP...\n")
immune <- RunUMAP(immune,
                             reduction = "harmony",
                             dims = 1:25,
                             verbose = FALSE)

cat("Clusters found:", length(levels(immune$seurat_clusters)), "\n")
print(table(immune$seurat_clusters))

# ── Step 5: Quick QC plots ────────────────────────────────────────────────────
p1 <- DimPlot(immune, group.by = "seurat_clusters",
              label = TRUE, label.size = 3.5, repel = TRUE) +
  ggtitle("All immune — clusters") + theme_pub() +
  theme(legend.position = "none")

p2 <- DimPlot(immune, group.by = "species",
              cols = species_colors) +
  ggtitle("Species") + theme_pub()

p3 <- DimPlot(immune,
              group.by = "immune_subtype",   # from prior annotation
              label = TRUE, label.size = 3, repel = TRUE) +
  ggtitle("Subtype (from primate/rodent annotation)") + theme_pub() +
  theme(legend.text = element_text(size = 8))

pdf(file.path(outdir, "05_umap_integrated_all_immune.pdf"), width = 20, height = 6)
print(p1 | p2 | p3)
dev.off()

cat("UMAP saved.\n")




# =============================================================================
# DotPlot + Cell Proportion Barplot — immune
# Groups: immune_subtype | Species
# =============================================================================

library(ggplot2)
library(dplyr)
library(patchwork)

# ── Marker genes from PanglaoDB / Tabula Muris / Tabula Sapiens ──────────────
# Ordered by lineage for a clean DotPlot layout
dotplot_genes_ordered <- c(
  # Pan leukocyte
  "PTPRC",
  # T cells
  "CD3E", "CD3D", "CD4", "TCF7", "CCR7", "SELL",
  "CD8A", "CD8B",
  "FOXP3", "IL2RA", "CTLA4",
  "GZMB", "GZMK", "PRF1",
  "PDCD1", "LAG3", "HAVCR2",
  "TRGC1", "TRDV1", "KLRB1",
  "RORC", "CXCR5", "BCL6",
  # NK
  "NCAM1", "NKG7", "GNLY", "KLRD1", "KLRF1", "TYROBP",
  # B cells
  "CD19", "MS4A1", "CD79A", "CD79B",
  "IGHD", "IGHM", "TCL1A",
  "CD27", "FCER2",
  "MZB1", "SDC1", "PRDM1", "XBP1", "IGHG1",
  # Monocytes
  "CD14", "LYZ", "S100A8", "S100A9", "VCAN",
  "FCGR3A", "CX3CR1",
  # Macrophages
  "CD68", "MRC1", "MARCO", "C1QA", "C1QB", "FOLR2",
  "TREM2", "APOE", "TNF", "IL1B",
  # DCs
  "ITGAX", "HLA-DRA",
  "CLEC9A", "XCR1",
  "CD1C", "CLEC10A", "FCER1A",
  "LILRA4", "CLEC4C", "IL3RA", "SIGLEC6",
  "CCR7", "CCL17", "CD83",
  # Mast
  "TPSAB1", "TPSB2", "KIT", "CPA3", "HDC", "GATA2", "MS4A2", "FCER1G",
  # Neutrophils
  "CXCR1", "CXCR2", "CSF3R", "FCGR3B", "MPO", "ELANE"
)

# Keep only genes present in the object
genes_ok <- intersect(dotplot_genes_ordered, rownames(immune))
cat("Genes found in object:", length(genes_ok), "of", length(dotplot_genes_ordered), "\n")

# ── Cell type order (lineage-sorted for DotPlot Y axis) ──────────────────────
subtype_order <- c(
  "CD4+ T cell", "Naive T cell", "Treg", "Th17", "Tfh",
  "CD8+ T cell", "Effector T cell", "Exhausted T cell",
  "gdT cell", "NKT cell",
  "NK cell",
  "B cell", "Naive B cell", "Memory B cell", "Plasma cell",
  "Classical monocyte", "Non-classical monocyte",
  "Macrophage", "Tissue-resident macrophage", "LAM macrophage", "M1 macrophage",
  "cDC1", "cDC2", "pDC", "Migratory DC",
  "Neutrophil",
  "Mast cell",
  "Unannotated"
)

# Only keep levels that exist in your object
subtype_order <- intersect(subtype_order, unique(immune$immune_subtype))

# Set factor order
immune$immune_subtype <- factor(immune$immune_subtype,
                                           levels = subtype_order)

# ── Color palette ─────────────────────────────────────────────────────────────
subtype_palette <- c(
  "CD4+ T cell"                = "#4DAF4A",
  "Naive T cell"               = "#80B1D3",
  "Treg"                       = "#A65628",
  "Th17"                       = "#FCCDE5",
  "Tfh"                        = "#D9D9D9",
  "CD8+ T cell"                = "#377EB8",
  "Effector T cell"            = "#1D9E75",
  "Exhausted T cell"           = "#BEBADA",
  "gdT cell"                   = "#FDB462",
  "NKT cell"                   = "#B3DE69",
  "NK cell"                    = "#BC80BD",
  "B cell"                     = "#E41A1C",
  "Naive B cell"               = "#FB8072",
  "Memory B cell"              = "#FFED6F",
  "Plasma cell"                = "#FF7F00",
  "Classical monocyte"         = "#984EA3",
  "Non-classical monocyte"     = "#CAB2D6",
  "Macrophage"                 = "#6A3D9A",
  "Tissue-resident macrophage" = "#33A02C",
  "LAM macrophage"             = "#B2DF8A",
  "M1 macrophage"              = "#A6CEE3",
  "cDC1"                       = "#1F78B4",
  "cDC2"                       = "#FDBF6F",
  "pDC"                        = "#FF6666",
  "Migratory DC"               = "#E31A1C",
  "Neutrophil"                 = "#B15928",
  "Mast cell"                  = "#00BFC4",
  "Unannotated"                = "#CCCCCC"
)

species_colors <- c(
  Human    = "#1D9E75",
  Marmoset = "#377EB8",
  Mouse    = "#E41A1C",
  Rat      = "#984EA3"
)

# =============================================================================
# PLOT 1 — DotPlot: immune_subtype (rows) × marker genes (cols)
# =============================================================================
Idents(immune) <- "immune_subtype"

p_dot <- suppressWarnings(
  DotPlot(
    immune,
    features = genes_ok,
    cols     = c("#313695", "white", "#a50026"),
    col.min  = -2,
    col.max  =  2,
    dot.scale = 5,
    scale    = TRUE
  ) +
    RotatedAxis() +
    scale_color_gradient2(low = "#313695", mid = "white", high = "#a50026",
                          midpoint = 0, name = "Avg\nexpression") +
    labs(x = NULL, y = NULL,
         title = "Immune subtype marker genes",
         size  = "% expressed") +
    theme_pub(base_size = 9) +
    theme(
      axis.text.x      = element_text(angle = 45, hjust = 1, size = 7,
                                      face = "italic"),
      axis.text.y      = element_text(size = 8),
      legend.key.size  = unit(0.35, "cm"),
      legend.text      = element_text(size = 7),
      legend.title     = element_text(size = 8),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
      plot.title       = element_text(face = "bold", hjust = 0.5)
    )
)

pdf(file.path(outdir, "06_dotplot_immune_subtypes.pdf"), width = 22, height = 8)
print(p_dot)
dev.off()

cat("DotPlot saved.\n")

# =============================================================================
# PLOT 2 — Cell proportion barplot: species (X) × immune_subtype (fill)
#           with absolute cell counts labeled
# =============================================================================

# ── Build count + proportion table ───────────────────────────────────────────
prop_df <- immune@meta.data %>%
  filter(!is.na(immune_subtype)) %>%
  group_by(species, immune_subtype) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(species) %>%
  mutate(
    total      = sum(n),
    proportion = n / total * 100
  ) %>%
  ungroup() %>%
  mutate(
    immune_subtype = factor(immune_subtype, levels = subtype_order),
    species        = factor(species, levels = c("Human", "Marmoset", "Mouse", "Rat"))
  )

# ── Species total labels (top of each bar) ───────────────────────────────────
totals_df <- prop_df %>%
  distinct(species, total)

# ── Stacked barplot (proportion) with count labels inside segments ────────────
# Label segments only if they are large enough to be readable (>= 3%)
label_df <- prop_df %>%
  filter(proportion >= 3) %>%
  group_by(species) %>%
  mutate(label_y = cumsum(proportion) - proportion / 2) %>%
  ungroup()

p_bar <- ggplot(prop_df,
                aes(x = species, y = proportion,
                    fill = immune_subtype)) +
  geom_bar(stat = "identity", width = 0.75, color = "white", linewidth = 0.2) +
  geom_text(data  = label_df,
            aes(y = label_y, label = n),
            size  = 2.5, color = "white", fontface = "bold") +
  geom_text(data  = totals_df,
            aes(x = species, y = 102, label = paste0("n = ", scales::comma(total)),
                fill = NULL),
            size  = 3, vjust = 0, fontface = "bold", color = "black") +
  scale_fill_manual(values = subtype_palette, name = "Immune subtype") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.06)),
                     labels = function(x) paste0(x, "%")) +
  scale_x_discrete(labels = c("Human", "Marmoset", "Mouse", "Rat")) +
  labs(x = NULL, y = "Cell proportion (%)",
       title = "Immune subtype composition per species") +
  theme_pub(base_size = 11) +
  theme(
    legend.position  = "right",
    legend.text      = element_text(size = 8),
    legend.key.size  = unit(0.45, "cm"),
    axis.text.x      = element_text(size = 10, face = "bold"),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3)
  ) +
  guides(fill = guide_legend(ncol = 1, reverse = FALSE))

pdf(file.path(outdir, "07_barplot_immune_proportion_by_species.pdf"),
    width = 10, height = 7)
print(p_bar)
dev.off()

cat("Barplot saved.\n")

# ── Also save the count table as CSV ─────────────────────────────────────────
write.csv(prop_df, file.path(outdir, "immune_subtype_counts_by_species.csv"),
          row.names = FALSE)
cat("Count table saved.\n")



# ── Genes not found ───────────────────────────────────────────────────────────
missing_genes <- setdiff(dotplot_genes_ordered, rownames(immune))
cat("Missing genes (", length(missing_genes), "):\n")
print(missing_genes)

# ── Reason 1: gene exists but under a different case/symbol ──────────────────
# Check case-insensitive match
rownames_upper <- toupper(rownames(immune))
missing_upper  <- toupper(missing_genes)

case_matches <- missing_genes[missing_upper %in% rownames_upper]
if (length(case_matches) > 0) {
  cat("\n── Case mismatch (gene exists under different capitalization) ──\n")
  for (g in case_matches) {
    actual <- rownames(immune)[toupper(rownames(immune)) == toupper(g)]
    cat(sprintf("  %s  →  found as: %s\n", g, paste(actual, collapse = ", ")))
  }
} else {
  cat("\nNo case mismatches found.\n")
}

# ── Reason 2: gene was filtered out during QC (low expression) ───────────────
# Check the raw counts slot of the ORIGINAL immune object (before HVG filtering)
# Note: ScaleData only keeps HVGs in scale.data, but counts/data keep all genes
cat("\n── Checking if missing genes are in raw counts (immune object) ──\n")
if ("counts" %in% Layers(immune[["RNA"]])) {
  raw_genes <- rownames(LayerData(immune, layer = "counts"))
} else {
  # fallback for joined layers
  raw_genes <- rownames(immune)
}
in_raw    <- intersect(missing_genes, raw_genes)
not_in_raw <- setdiff(missing_genes, raw_genes)

cat("Present in raw counts but not in scaled/HVG slot:", length(in_raw), "\n")
print(in_raw)
cat("\nNot in raw counts at all (filtered in QC or absent in all species):", length(not_in_raw), "\n")
print(not_in_raw)

# ── Reason 3: gene is species-specific (check per species) ───────────────────
cat("\n── Checking per-species presence of missing genes ──\n")
species_list <- c("Human", "Marmoset", "Mouse", "Rat")

presence_mat <- sapply(species_list, function(sp) {
  cells_sp <- colnames(immune)[immune$species == sp]
  sub_obj   <- subset(immune, cells = cells_sp)
  not_in_raw %>%
    sapply(function(g) g %in% rownames(sub_obj))
})

rownames(presence_mat) <- not_in_raw
cat("\nGene presence per species (TRUE = present):\n")
print(presence_mat)

# ── Reason 4: synonym / alias — check a few known aliases ────────────────────
cat("\n── Common alias check ──\n")
aliases <- c(
  "HAVCR2"  = "TIM3",
  "MS4A1"   = "CD20",
  "NCAM1"   = "CD56",
  "FCGR3A"  = "CD16A",
  "PDCD1"   = "PD1",
  "ITGAX"   = "CD11C",
  "SDC1"    = "CD138",
  "IL2RA"   = "CD25",
  "FCER2"   = "CD23",
  "SELL"    = "CD62L"
)

for (canonical in names(aliases)) {
  alias <- aliases[canonical]
  if (canonical %in% missing_genes) {
    found_alias <- alias %in% rownames(immune)
    cat(sprintf("  %s (alias: %s) → alias found in object: %s\n",
                canonical, alias, found_alias))
  }
}




# =============================================================================
# DotPlot per species — immune subtypes × marker genes
# Uses immune with immune_subtype metadata
# Handles mixed gene namespaces across species
# =============================================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

outdir <- "immune_reintegration"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

theme_pub <- function(base_size = 11) {
  theme_classic(base_size = base_size, base_family = "") +
    theme(
      axis.line         = element_line(linewidth = 0.4),
      axis.ticks        = element_line(linewidth = 0.3),
      axis.text         = element_text(color = "black"),
      legend.background = element_blank(),
      strip.background  = element_blank(),
      strip.text        = element_text(face = "bold"),
      plot.title        = element_text(face = "bold", hjust = 0.5)
    )
}

# =============================================================================
# STEP 1 — Marker genes tailored to your 13 subtypes
# These are the canonical markers from PanglaoDB/Tabula Muris/Tabula Sapiens
# mapped to your actual subtype labels
# =============================================================================

subtype_markers <- list(
  "Naive/CM T cells"              = c("CD3E", "CD3D", "CD4", "TCF7",
                                      "CCR7", "SELL", "IL7R", "LEF1"),
  "NK/Cytotoxic T cells"          = c("CD8A", "NKG7", "GZMK", "PRF1",
                                      "NCAM1", "TYROBP", "KLRD1", "GNLY"),
  "gdT cells"                     = c("CD3E", "TRGC1", "TRDV1", "KLRB1"),
  "B cells"                       = c("CD19", "MS4A1", "CD79A", "CD79B",
                                      "IGHD", "TCL1A"),
  "Plasma cells"                  = c("MZB1", "SDC1", "XBP1", "PRDM1",
                                      "IGHG1", "JCHAIN"),
  "Proliferating Plasma cells"    = c("MZB1", "XBP1", "TOP2A", "MKI67",
                                      "PCNA", "STMN1"),
  "Scavenger Macrophages"         = c("CD68", "MSR1", "MRC1", "APOE",
                                      "TREM2", "C1QA", "C1QB", "FOLR2"),
  "Tissue-resident Macrophages"   = c("CD68", "MARCO", "C1QA", "C1QB",
                                      "FOLR2", "LYVE1", "TIMD4"),
  "AXL+ Dendritic cells"          = c("AXL", "SIGLEC6", "LILRA4", "CLEC4C",
                                      "CLEC9A", "XCR1", "CD1C", "ITGAX"),
  "Classical Neutrophils"         = c("S100A8", "S100A9", "CSF3R", "FCGR3B",
                                      "MPO", "ELANE", "CXCR2"),
  "Inflammatory Neutrophils"      = c("S100A8", "S100A9", "IL1B", "TNF",
                                      "CXCL8", "CCL3", "NLRP3"),
  "Mast Cell Progenitors"         = c("KIT", "GATA2", "FCER1G", "MS4A2",
                                      "CPA3", "HDC"),
  "Mature Mast Cells (MCTC)"      = c("TPSAB1", "TPSB2", "KIT", "CPA3",
                                      "MS4A2", "FCER1G", "HPGDS")
)

# Ordered subtypes matching your barplot
subtype_order <- c(
  "Naive/CM T cells",
  "NK/Cytotoxic T cells",
  "\u03b3\u03b4 T cells",           # γδ T cells — unicode to avoid encoding issues
  "B cells",
  "Plasma cells",
  "Proliferating Plasma cells",
  "Scavenger Macrophages",
  "Tissue-resident Macrophages",
  "AXL+ Dendritic cells",
  "Classical Neutrophils",
  "Inflammatory Neutrophils",
  "Mast Cell Progenitors",
  "Mature Mast Cells (MCTC)"
)

# Keep only levels present in object
subtype_order_present <- intersect(subtype_order, unique(immune$immune_subtype))

# =============================================================================
# STEP 2 — Check which genes are available per species
# =============================================================================
cat("── Checking gene availability per species ──\n")

species_list  <- c("Human", "Marmoset", "Mouse", "Rat")
all_markers   <- unique(unlist(subtype_markers))

gene_presence <- sapply(species_list, function(sp) {
  cells_sp <- colnames(immune)[immune$species == sp]
  rn       <- rownames(immune[, cells_sp])
  all_markers %in% rn
})
rownames(gene_presence) <- all_markers

cat("\nGene presence across species:\n")
print(gene_presence)

# Genes present in at least 2 species
genes_keep <- all_markers[rowSums(gene_presence) >= 2]
cat("\nGenes present in ≥2 species:", length(genes_keep), "\n")

# =============================================================================
# STEP 3 — Build per-species DotPlots
# =============================================================================

# Set factor levels
immune$immune_subtype <- factor(
  immune$immune_subtype,
  levels = subtype_order_present
)

# Ordered gene list for DotPlot (preserve lineage order, filter to available)
genes_ordered <- unique(unlist(subtype_markers))
genes_ordered <- intersect(genes_ordered, genes_keep)

# Subtype color palette
subtype_palette <- c(
  "Naive/CM T cells"              = "#80B1D3",
  "NK/Cytotoxic T cells"          = "#4DAF4A",
  "\u03b3\u03b4 T cells"          = "#FDB462",
  "B cells"                       = "#E41A1C",
  "Plasma cells"                  = "#FF7F00",
  "Proliferating Plasma cells"    = "#FFED6F",
  "Scavenger Macrophages"         = "#6A3D9A",
  "Tissue-resident Macrophages"   = "#33A02C",
  "AXL+ Dendritic cells"          = "#1F78B4",
  "Classical Neutrophils"         = "#B15928",
  "Inflammatory Neutrophils"      = "#A6CEE3",
  "Mast Cell Progenitors"         = "#00BFC4",
  "Mature Mast Cells (MCTC)"      = "#984EA3"
)

# ── Function: DotPlot for one species ────────────────────────────────────────
make_species_dotplot <- function(seu, species_name, genes, subtype_ord,
                                 palette, base_size = 9) {
  
  # Subset to species
  cells_sp <- colnames(seu)[seu$species == species_name]
  sub_obj  <- subset(seu, cells = cells_sp)
  
  # Filter to genes present in this species
  genes_sp <- intersect(genes, rownames(sub_obj))
  cat(sprintf("  [%s] %d / %d genes available\n",
              species_name, length(genes_sp), length(genes)))
  
  if (length(genes_sp) == 0) {
    cat(sprintf("  [%s] No genes found — skipping\n", species_name))
    return(NULL)
  }
  
  # Keep subtypes present in this species
  subtypes_sp <- intersect(subtype_ord,
                           unique(sub_obj$immune_subtype[!is.na(sub_obj$immune_subtype)]))
  sub_obj$immune_subtype <- factor(sub_obj$immune_subtype, levels = subtypes_sp)
  
  Idents(sub_obj) <- "immune_subtype"
  
  p <- suppressWarnings(
    DotPlot(
      sub_obj,
      features  = genes_sp,
      dot.scale = 5,
      scale     = TRUE,
      col.min   = -2,
      col.max   =  2
    ) +
      RotatedAxis() +
      scale_color_gradient2(
        low      = "#313695",
        mid      = "white",
        high     = "#a50026",
        midpoint = 0,
        name     = "Avg\nexp"
      ) +
      labs(
        x     = NULL,
        y     = NULL,
        title = species_name,
        size  = "% exp"
      ) +
      theme_pub(base_size = base_size) +
      theme(
        axis.text.x      = element_text(angle = 45, hjust = 1,
                                        size = 6.5, face = "italic"),
        axis.text.y      = element_text(size = 7.5),
        legend.key.size  = unit(0.3, "cm"),
        legend.text      = element_text(size = 6.5),
        legend.title     = element_text(size = 7),
        panel.grid.major = element_line(color = "grey92", linewidth = 0.25),
        plot.title       = element_text(face = "bold", hjust = 0.5,
                                        size = 10, color = "black")
      )
  )
  return(p)
}

# ── Generate one panel per species ───────────────────────────────────────────
cat("\n── Generating per-species DotPlots ──\n")

plots <- lapply(species_list, function(sp) {
  make_species_dotplot(
    seu         = immune,
    species_name = sp,
    genes        = genes_ordered,
    subtype_ord  = subtype_order_present,
    palette      = subtype_palette
  )
})
names(plots) <- species_list

# Remove NULLs
plots <- Filter(Negate(is.null), plots)

# =============================================================================
# STEP 4 — Save individual PDFs per species
# =============================================================================
cat("\n── Saving per-species DotPlot PDFs ──\n")

for (sp in names(plots)) {
  fname <- file.path(outdir,
                     paste0("09_dotplot_", tolower(sp), "_immune_subtypes.pdf"))
  pdf(fname, width = 16, height = 6)
  print(plots[[sp]])
  dev.off()
  cat(sprintf("  Saved: %s\n", fname))
}

# =============================================================================
# STEP 5 — Combined 4-panel figure (stacked vertically)
# =============================================================================
cat("\n── Saving combined 4-panel DotPlot ──\n")

combined <- wrap_plots(plots, ncol = 1) +
  plot_annotation(
    title   = "Immune subtype marker expression per species",
    theme   = theme(plot.title = element_text(face = "bold",
                                              hjust = 0.5, size = 13))
  )

pdf(file.path(outdir, "09_dotplot_all_species_immune_subtypes.pdf"),
    width = 18, height = 22)
print(combined)
dev.off()

cat("Combined DotPlot saved.\n")

# =============================================================================
# STEP 6 — Fix γδ T cells label in barplot (bonus fix)
# Replace the unicode character issue by using "gdT cells" as safe fallback
# or export with cairo_pdf which handles unicode properly
# =============================================================================
cat("\n── Tip: to fix Greek character in PDF, use cairo_pdf() instead of pdf() ──\n")
cat("Example:\n")
cat('  cairo_pdf(file.path(outdir, "07_barplot_fixed.pdf"), width = 9, height = 7)\n')
cat("  print(p_bar)\n")
cat("  dev.off()\n")


# =============================================================================
# DotPlot with gene-subtype annotation bar on top
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# ── Build gene → subtype annotation table ────────────────────────────────────
gene_annotation <- stack(subtype_markers) %>%
  rename(gene = values, subtype = ind) %>%
  mutate(subtype = as.character(subtype)) %>%
  # If a gene appears in multiple subtypes, keep first occurrence (lineage order)
  group_by(gene) %>%
  slice(1) %>%
  ungroup() %>%
  filter(gene %in% genes_ordered)

# Match gene order to dotplot
gene_annotation$gene <- factor(gene_annotation$gene, levels = genes_ordered)
gene_annotation$subtype <- factor(gene_annotation$subtype,
                                  levels = names(subtype_markers))

# ── Color palette for gene annotation bar ────────────────────────────────────
# One color per subtype label on the x-axis annotation
annotation_palette <- c(
  "Naive/CM T cells"              = "#80B1D3",
  "NK/Cytotoxic T cells"          = "#4DAF4A",
  "gdT cells"                     = "#FDB462",
  "B cells"                       = "#E41A1C",
  "Plasma cells"                  = "#FF7F00",
  "Proliferating Plasma cells"    = "#FFED6F",
  "Scavenger Macrophages"         = "#6A3D9A",
  "Tissue-resident Macrophages"   = "#33A02C",
  "AXL+ Dendritic cells"          = "#1F78B4",
  "Classical Neutrophils"         = "#B15928",
  "Inflammatory Neutrophils"      = "#A6CEE3",
  "Mast Cell Progenitors"         = "#00BFC4",
  "Mature Mast Cells (MCTC)"      = "#984EA3"
)

# ── Function: extract DotPlot data + add annotation bar ──────────────────────
make_annotated_dotplot <- function(seu, species_name, genes,
                                   subtype_ord, ann_palette,
                                   gene_ann, base_size = 9) {
  
  # Subset to species
  cells_sp <- colnames(seu)[seu$species == species_name]
  sub_obj  <- subset(seu, cells = cells_sp)
  genes_sp <- intersect(genes, rownames(sub_obj))
  
  subtypes_sp <- intersect(subtype_ord,
                           unique(sub_obj$immune_subtype[!is.na(sub_obj$immune_subtype)]))
  sub_obj$immune_subtype <- factor(sub_obj$immune_subtype, levels = subtypes_sp)
  Idents(sub_obj) <- "immune_subtype"
  
  cat(sprintf("  [%s] %d genes | %d subtypes\n",
              species_name, length(genes_sp), length(subtypes_sp)))
  
  # ── Extract DotPlot data ──────────────────────────────────────────────────
  dot_data <- suppressWarnings(
    DotPlot(sub_obj, features = genes_sp, scale = TRUE,
            col.min = -2, col.max = 2)$data
  )
  
  # Rename columns for clarity
  dot_data <- dot_data %>%
    rename(gene = features.plot, subtype_y = id) %>%
    mutate(
      gene      = factor(gene, levels = genes_sp),
      subtype_y = factor(subtype_y, levels = subtypes_sp)
    )
  
  # ── Main DotPlot panel ────────────────────────────────────────────────────
  p_main <- ggplot(dot_data,
                   aes(x = gene, y = subtype_y,
                       size = pct.exp, color = avg.exp.scaled)) +
    geom_point() +
    scale_color_gradient2(low = "#313695", mid = "white", high = "#a50026",
                          midpoint = 0, name = "Avg\nexp",
                          limits = c(-2, 2), oob = scales::squish) +
    scale_size_continuous(name = "% exp", range = c(0.2, 5),
                          breaks = c(10, 25, 50, 75, 100)) +
    scale_x_discrete(position = "bottom") +
    labs(x = NULL, y = NULL, title = species_name) +
    theme_pub(base_size = base_size) +
    theme(
      axis.text.x      = element_text(angle = 45, hjust = 1,
                                      size = 6.5, face = "italic"),
      axis.text.y      = element_text(size = 8),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.25),
      legend.key.size  = unit(0.3, "cm"),
      legend.text      = element_text(size = 6.5),
      legend.title     = element_text(size = 7),
      plot.title       = element_text(face = "bold", hjust = 0.5, size = 10),
      plot.margin      = margin(t = 2, r = 5, b = 2, l = 5)
    )
  
  # ── Annotation bar (top strip) ────────────────────────────────────────────
  ann_data <- gene_ann %>%
    filter(gene %in% genes_sp) %>%
    mutate(gene = factor(gene, levels = genes_sp),
           y    = 1)
  
  p_ann <- ggplot(ann_data, aes(x = gene, y = y, fill = subtype)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_manual(values  = ann_palette,
                      name    = "Marker\nsubtype",
                      drop    = FALSE) +
    scale_x_discrete(position = "top") +
    scale_y_continuous(expand = c(0, 0)) +
    theme_void(base_family = "") +
    theme(
      axis.text.x      = element_blank(),
      legend.position  = "none",          # legend shown in main plot
      plot.margin      = margin(t = 2, r = 5, b = 0, l = 5)
    )
  
  # ── Annotation label strip (subtype names under tiles) ───────────────────
  # Find midpoint position of each subtype block on x-axis
  label_data <- ann_data %>%
    mutate(gene_num = as.integer(gene)) %>%
    group_by(subtype) %>%
    summarise(
      xmid  = mean(gene_num),
      xmin  = min(gene_num) - 0.5,
      xmax  = max(gene_num) + 0.5,
      .groups = "drop"
    ) %>%
    mutate(y = 1)
  
  p_label <- ggplot(label_data,
                    aes(x = xmid, y = y,
                        label = subtype, color = subtype)) +
    geom_text(size = 2.2, fontface = "bold", hjust = 0.5, vjust = 0.5,
              show.legend = FALSE) +
    geom_segment(aes(x = xmin, xend = xmax, y = 0.45, yend = 0.45,
                     color = subtype), linewidth = 0.5, show.legend = FALSE) +
    scale_color_manual(values = ann_palette) +
    scale_x_continuous(limits = c(0.5, length(genes_sp) + 0.5),
                       expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1.5), expand = c(0, 0)) +
    theme_void(base_family = "") +
    theme(plot.margin = margin(t = 0, r = 5, b = 1, l = 5))
  
  # ── Stack: annotation bar / label strip / main dotplot ───────────────────
  p_ann / p_label / p_main +
    plot_layout(heights = c(0.06, 0.1, 1))
}

# ── Generate plots ────────────────────────────────────────────────────────────
cat("\n── Generating annotated DotPlots per species ──\n")

plots_ann <- lapply(species_list, function(sp) {
  make_annotated_dotplot(
    seu          = immune,
    species_name = sp,
    genes        = genes_ordered,
    subtype_ord  = subtype_order_present,
    ann_palette  = annotation_palette,
    gene_ann     = gene_annotation
  )
})
names(plots_ann) <- species_list
plots_ann <- Filter(Negate(is.null), plots_ann)

# ── Save individual PDFs ──────────────────────────────────────────────────────
for (sp in names(plots_ann)) {
  fname <- file.path(getwd(),
                     paste0("10_dotplot_annotated_", tolower(sp), ".pdf"))
  cairo_pdf(fname, width = 16, height = 7)
  print(plots_ann[[sp]])
  dev.off()
  cat(sprintf("  Saved: %s\n", fname))
}

# ── Combined 4-panel ──────────────────────────────────────────────────────────
cairo_pdf(file.path(getwd(), "10_dotplot_annotated_all_species.pdf"),
          width = 18, height = 26)
print(
  wrap_plots(plots_ann, ncol = 1) +
    plot_annotation(
      title = "Immune subtype marker expression per species",
      theme = theme(plot.title = element_text(face = "bold",
                                              hjust = 0.5, size = 13))
    )
)
dev.off()
cat("Done.\n")

### PAREI AQUI
library(cowplot)

# ── UMAP 1: all cells, coloured by epi_subtype ────────────────────────────────
epi_clean$epi_subtype <- factor(epi_clean$epi_subtype, levels = subtype_order_present)

p_all <- DimPlot(epi_clean, group.by = "epi_subtype",
                 cols = subtype_palette,
                 label = TRUE, label.size = 3.5, repel = TRUE, pt.size = 0.3) +
  ggtitle("Epithelial subtypes — all species") +
  theme_pub() +
  theme(legend.text = element_text(size = 9),
        legend.key.size = unit(0.45, "cm"))

cairo_pdf(file.path(outdir, "08_umap_epi_subtypes_all.pdf"), width = 10, height = 8)
print(p_all)
dev.off()

# ── UMAP 2: split by species ───────────────────────────────────────────────────
panels <- lapply(c("Human", "Marmoset", "Mouse", "Rat"), function(sp) {
  n_cells <- sum(epi_clean$species == sp, na.rm = TRUE)
  DimPlot(epi_clean,
          group.by  = "epi_subtype",
          cells     = colnames(epi_clean)[epi_clean$species == sp],
          cols      = subtype_palette,
          label     = FALSE, pt.size = 0.5) +
    ggtitle(paste0(sp, "\n(n = ", scales::comma(n_cells), ")")) +
    theme_pub(base_size = 9) +
    theme(legend.position = "none",
          axis.text  = element_blank(),
          axis.ticks = element_blank(),
          axis.line  = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5, size = 10))
})

# Shared legend from the combined plot
legend_plot <- DimPlot(epi_clean, group.by = "epi_subtype", cols = subtype_palette) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3))) +
  theme(legend.text = element_text(size = 9),
        legend.key.size = unit(0.45, "cm"))
shared_legend <- cowplot::get_legend(legend_plot)

cairo_pdf(file.path(outdir, "09_umap_epi_subtypes_split_species.pdf"), width = 22, height = 6)
print((panels[[1]] | panels[[2]] | panels[[3]] | panels[[4]]) |
        wrap_elements(shared_legend))
dev.off()

cat("Saved: 08_umap_epi_subtypes_all.pdf\n")
cat("Saved: 09_umap_epi_subtypes_split_species.pdf\n")






# ── Step 6: Save ──────────────────────────────────────────────────────────────
saveRDS(immune, file.path(outdir, "immune_all4.rds"))
cat("Saved: immune_all4.rds\n")
