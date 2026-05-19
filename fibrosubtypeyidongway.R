###fibro sub typing ###

setwd("/r_workspace/cristal_data/cross_species/subtypes_cross_species/fibroblast")

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(scales)
})

set.seed(42)

# Quick check
cat("Object name: fibroblast\n")
cat("Total cells:", ncol(fibroblast), "\n")
cat("Species:\n")
print(table(fibroblast$species))
cat("\nAssays:", Assays(fibroblast), "\n")
cat("Layers:", Layers(fibroblast[["RNA"]]), "\n")
cat("\nMetadata columns:\n")
print(colnames(fibroblast@meta.data))


# =============================================================================
# Fibroblast Compartment Re-integration & Subtype Annotation
# Input:    `fibroblast` Seurat v5 object (4 species, 22,697 cells)
# Strategy: Human + Marmoset → Harmony | Mouse + Rat → Harmony → annotate
# Markers:  PanglaoDB + Tabula Sapiens + Tabula Muris (fibroblast/stromal)
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)
  library(cowplot)
})

set.seed(42)

# ─── Settings ─────────────────────────────────────────────────────────────────
outdir <- "fibroblast_reintegration"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

RESOLUTION_PRIMATE <- 0.5
RESOLUTION_RODENT  <- 0.5
N_PCS              <- 30
N_HARMONY_DIMS     <- 25

species_colors <- c(
  Human    = "#1D9E75",
  Marmoset = "#377EB8",
  Mouse    = "#E41A1C",
  Rat      = "#984EA3"
)

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
# STEP 0 — Confirm object
# =============================================================================
cat("── Step 0: Confirming fibroblast object ──\n")
cat("Total cells:", ncol(fibroblast), "\n")
cat("Species breakdown:\n")
print(table(fibroblast$species))
cat("Metadata columns:\n")
print(colnames(fibroblast@meta.data))

# Drop leftover assays to keep object clean
DefaultAssay(fibroblast) <- "RNA"
if ("SCT" %in% Assays(fibroblast))        fibroblast[["SCT"]]        <- NULL
if ("integrated" %in% Assays(fibroblast)) fibroblast[["integrated"]] <- NULL
fibroblast[["RNA"]] <- JoinLayers(fibroblast[["RNA"]])
cat("Assays after cleanup:", Assays(fibroblast), "\n")

# =============================================================================
# STEP 1 — Split into primate and rodent
# =============================================================================
cat("\n── Step 1: Splitting into primate and rodent groups ──\n")

fib_primate <- subset(fibroblast, subset = species %in% c("Human", "Marmoset"))
fib_rodent  <- subset(fibroblast, subset = species %in% c("Mouse", "Rat"))

cat("Primate — Human:", sum(fib_primate$species == "Human"),
    "| Marmoset:", sum(fib_primate$species == "Marmoset"),
    "| Total:", ncol(fib_primate), "\n")
cat("Rodent  — Mouse:", sum(fib_rodent$species == "Mouse"),
    "| Rat:", sum(fib_rodent$species == "Rat"),
    "| Total:", ncol(fib_rodent), "\n")

# =============================================================================
# STEP 2 — Pre-processing helper (NormalizeData → PCA)
# =============================================================================
preprocess_seurat <- function(seu, group_name, n_pcs = N_PCS) {
  cat(sprintf("\n  [%s] Normalizing...\n", group_name))
  seu <- NormalizeData(seu, verbose = FALSE)
  
  cat(sprintf("  [%s] Finding variable features (split by species)...\n", group_name))
  seu[["RNA"]] <- JoinLayers(seu[["RNA"]])
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
fib_primate <- preprocess_seurat(fib_primate, "Primate")
fib_primate <- run_harmony(fib_primate, "Primate", RESOLUTION_PRIMATE)
saveRDS(fib_primate, file.path(outdir, "fib_primate.rds"))

fib_rodent  <- preprocess_seurat(fib_rodent,  "Rodent")
fib_rodent  <- run_harmony(fib_rodent,  "Rodent",  RESOLUTION_RODENT)
saveRDS(fib_rodent, file.path(outdir, "fib_rodent.rds"))

cat("\n  Both objects integrated and saved.\n")

# =============================================================================
# STEP 4 — Fibroblast subtype marker genes
#   Sources: PanglaoDB (Franzén et al. 2019)
#            Tabula Sapiens (2022 Science) — tongue/oral connective tissue
#            Tabula Muris (2018 Nature)    — tongue fibroblasts
#            Oral/tongue CAF literature (Bhattacharya et al., Givel et al.)
# =============================================================================

subtype_markers <- list(
  
  # ── Universal fibroblast identity ────────────────────────────────────────
  # PanglaoDB: pan-fibroblast markers
  "Pan-fibroblast"            = c("DCN", "LUM", "COL1A1", "COL1A2",
                                  "COL3A1", "VIM", "PDGFRA", "FBN1"),
  
  # ── Resident / homeostatic fibroblasts ───────────────────────────────────
  # Tabula Sapiens: resting stromal fibroblasts; tongue lamina propria
  "Resident fibroblast"       = c("MFAP5", "MFAP4", "CFD", "GSN",
                                  "FBLN1", "FBLN2", "SFRP2", "SFRP4"),
  
  # ── Adventitial / CD34+ fibroblasts ──────────────────────────────────────
  # PanglaoDB + Tabula Sapiens: perivascular niche fibroblasts
  "Adventitial fibroblast"    = c("CD34", "PI16", "LY6C1", "CXCL12",
                                  "ADGRE5", "ENTPD2", "IL33", "DPP4"),
  
  # ── Myofibroblasts ────────────────────────────────────────────────────────
  # PanglaoDB + Tabula Muris: activated/contractile fibroblasts
  "Myofibroblast"             = c("ACTA2", "TAGLN", "MYH11", "CNN1",
                                  "TPM2", "MMP11", "POSTN", "LOXL2"),
  
  # ── Inflammatory CAFs ─────────────────────────────────────────────────────
  # Givel et al. 2018; Bhattacharya oral CAF literature
  "Inflammatory CAF"          = c("IL6", "IL8", "CXCL1", "CXCL2",
                                  "CXCL8", "CCL2", "LIF", "PTGES"),
  
  # ── Antigen-presenting CAFs ───────────────────────────────────────────────
  # Tabula Sapiens: MHC-II+ fibroblasts; oral immune-stromal niche
  "Antigen-presenting CAF"    = c("HLA-DRA", "HLA-DRB1", "HLA-DPA1",
                                  "CD74", "CIITA", "HLA-DQA1"),
  
  # ── ECM-remodeling fibroblasts ───────────────────────────────────────────
  # Tabula Muris + oral literature: matrix-producing subtype
  "ECM-remodeling fibroblast" = c("MMP1", "MMP3", "MMP13", "TIMP1",
                                  "TIMP3", "LOX", "LOXL1", "COL4A1"),
  
  # ── Lipofibroblasts ───────────────────────────────────────────────────────
  # Tabula Muris: lipid-storing fibroblasts; tongue submucosal layer
  "Lipofibroblast"            = c("PLIN2", "PLIN1", "ADIPOQ", "FABP4",
                                  "LPL", "PPARG", "CEBPA", "FASN"),
  
  # ── Pericytes ─────────────────────────────────────────────────────────────
  # PanglaoDB: mural cells; often co-isolated with fibroblasts
  "Pericyte"                  = c("RGS5", "PDGFRB", "NOTCH3", "MCAM",
                                  "CSPG4", "ABCC9", "KCNJ8", "GJA4"),
  
  # ── Smooth muscle cells ───────────────────────────────────────────────────
  # PanglaoDB: vascular smooth muscle; tongue vasculature
  "Smooth muscle"             = c("ACTA2", "MYH11", "CNN1", "SMTN",
                                  "MYLK", "ACTG2", "TAGLN", "DES"),
  
  # ── Proliferating fibroblasts ─────────────────────────────────────────────
  # Tabula Sapiens + Tabula Muris: cycling stromal cells
  "Proliferating"             = c("MKI67", "TOP2A", "PCNA", "STMN1",
                                  "CENPF", "UBE2C", "BIRC5", "CDC20"),
  
  # ── Neural crest-derived fibroblasts ─────────────────────────────────────
  # Tongue-specific; neural crest contribution to oral mesenchyme
  "Neural crest fibroblast"   = c("SOX10", "S100B", "TWIST1", "SNAI2",
                                  "FOXD3", "PAX3", "NES", "NGFR")
)

# ── Ordered gene list ─────────────────────────────────────────────────────────
genes_ordered <- unique(unlist(subtype_markers))

# ── Subtype order ─────────────────────────────────────────────────────────────
subtype_order <- c(
  "Pan-fibroblast",
  "Resident fibroblast",
  "Adventitial fibroblast",
  "Myofibroblast",
  "Inflammatory CAF",
  "Antigen-presenting CAF",
  "ECM-remodeling fibroblast",
  "Lipofibroblast",
  "Pericyte",
  "Smooth muscle",
  "Proliferating",
  "Neural crest fibroblast"
)

# ── Color palette ─────────────────────────────────────────────────────────────
subtype_palette <- c(
  "Pan-fibroblast"            = "#AAAAAA",
  "Resident fibroblast"       = "#E41A1C",
  "Adventitial fibroblast"    = "#FF7F00",
  "Myofibroblast"             = "#FFED6F",
  "Inflammatory CAF"          = "#4DAF4A",
  "Antigen-presenting CAF"    = "#80B1D3",
  "ECM-remodeling fibroblast" = "#377EB8",
  "Lipofibroblast"            = "#984EA3",
  "Pericyte"                  = "#A65628",
  "Smooth muscle"             = "#F781BF",
  "Proliferating"             = "#FFED6F",
  "Neural crest fibroblast"   = "#B2DF8A",
  "Contaminant"               = "#CCCCCC"
)

# =============================================================================
# STEP 5 — Module scores per subtype
# =============================================================================
cat("\n── Step 5: Module scores ──\n")

add_module_scores <- function(seu, markers, obj_name) {
  for (nm in names(markers)) {
    genes_present <- intersect(markers[[nm]], rownames(seu))
    if (length(genes_present) == 0) {
      cat(sprintf("  [%s] Skipping '%s' — no genes found\n", obj_name, nm))
      next
    }
    safe_nm <- gsub("[^A-Za-z0-9_]", "_", nm)
    seu <- AddModuleScore(seu,
                          features = list(genes_present),
                          name     = paste0("score_", safe_nm),
                          ctrl     = 50, seed = 42)
    colnames(seu@meta.data)[ncol(seu@meta.data)] <- paste0("score_", safe_nm)
  }
  return(seu)
}

fib_primate <- add_module_scores(fib_primate, subtype_markers, "Primate")
fib_rodent  <- add_module_scores(fib_rodent,  subtype_markers, "Rodent")

# =============================================================================
# STEP 6 — FindAllMarkers (save CSVs for inspection)
# =============================================================================
cat("\n── Step 6: FindAllMarkers ──\n")

Idents(fib_primate) <- "seurat_clusters"
markers_primate <- FindAllMarkers(fib_primate, only.pos = TRUE,
                                  min.pct = 0.25, logfc.threshold = 0.3,
                                  test.use = "wilcox", verbose = FALSE)
write.csv(markers_primate,
          file.path(outdir, "markers_fib_primate_clusters.csv"),
          row.names = FALSE)

Idents(fib_rodent) <- "seurat_clusters"
markers_rodent <- FindAllMarkers(fib_rodent, only.pos = TRUE,
                                 min.pct = 0.25, logfc.threshold = 0.3,
                                 test.use = "wilcox", verbose = FALSE)
write.csv(markers_rodent,
          file.path(outdir, "markers_fib_rodent_clusters.csv"),
          row.names = FALSE)

cat("  Marker CSVs saved.\n")

# =============================================================================
# STEP 7 — Exploratory plots (pre-annotation)
# =============================================================================
cat("\n── Step 7: Exploratory plots ──\n")

key_features <- c(
  "DCN", "LUM", "COL1A1",          # pan-fibroblast
  "CD34", "PI16", "CXCL12",        # adventitial
  "ACTA2", "TAGLN", "POSTN",       # myofibroblast
  "IL6", "CCL2", "CXCL1",          # inflammatory CAF
  "RGS5", "PDGFRB", "NOTCH3",      # pericyte
  "MKI67", "TOP2A",                # proliferating
  "PLIN2", "ADIPOQ",               # lipofibroblast
  "SOX10", "S100B"                 # neural crest
)

dotplot_genes <- c(
  "DCN", "LUM", "COL1A1", "COL1A2", "COL3A1", "VIM", "PDGFRA",
  "MFAP5", "FBLN1", "SFRP2",
  "CD34", "PI16", "CXCL12", "DPP4",
  "ACTA2", "TAGLN", "MYH11", "POSTN", "LOXL2",
  "IL6", "CCL2", "CXCL1", "LIF",
  "HLA-DRA", "CD74", "CIITA",
  "MMP1", "MMP3", "TIMP1", "LOX",
  "PLIN2", "ADIPOQ", "FABP4",
  "RGS5", "PDGFRB", "NOTCH3", "MCAM",
  "MYH11", "CNN1", "MYLK", "DES",
  "MKI67", "TOP2A", "PCNA",
  "SOX10", "S100B", "TWIST1", "NGFR"
)

make_plots <- function(seu, group_name, spc_colors, dot_genes, feat_genes) {
  
  p_clust <- DimPlot(seu, group.by = "seurat_clusters",
                     label = TRUE, label.size = 3.5, repel = TRUE) +
    ggtitle(paste(group_name, "— clusters")) + theme_pub() +
    theme(legend.position = "none")
  p_spc <- DimPlot(seu, group.by = "species", cols = spc_colors) +
    ggtitle("Species") + theme_pub()
  
  cairo_pdf(file.path(outdir,
                      paste0("01_umap_", tolower(gsub(" ", "_", group_name)),
                             "_exploratory.pdf")),
            width = 14, height = 6)
  print(p_clust | p_spc)
  dev.off()
  
  genes_ok <- intersect(dot_genes, rownames(seu))
  Idents(seu) <- "seurat_clusters"
  p_dot <- suppressWarnings(
    DotPlot(seu, features = genes_ok, dot.scale = 5) +
      RotatedAxis() +
      scale_color_gradient2(low = "#313695", mid = "white", high = "#a50026",
                            midpoint = 0) +
      ggtitle(paste(group_name, "— fibroblast subtype markers")) +
      theme_pub(base_size = 9) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7,
                                       face = "italic"),
            legend.key.size = unit(0.35, "cm"))
  )
  
  cairo_pdf(file.path(outdir,
                      paste0("02_dotplot_", tolower(gsub(" ", "_", group_name)),
                             ".pdf")),
            width = 20, height = 7)
  print(p_dot)
  dev.off()
  
  feats_ok <- intersect(feat_genes, rownames(seu))
  p_feat <- FeaturePlot(seu, features = feats_ok,
                        order = TRUE, ncol = 4,
                        cols = c("lightgrey", "#a50026"), pt.size = 0.3) &
    theme_pub(base_size = 9)
  
  cairo_pdf(file.path(outdir,
                      paste0("03_featureplot_", tolower(gsub(" ", "_", group_name)),
                             ".pdf")),
            width = 18, height = 12)
  print(p_feat)
  dev.off()
  
  cat(sprintf("  [%s] Exploratory plots saved.\n", group_name))
}

make_plots(fib_primate, "Primate", species_colors, dotplot_genes, key_features)
make_plots(fib_rodent,  "Rodent",  species_colors, dotplot_genes, key_features)

# =============================================================================
# STEP 8 — ANNOTATION
# ─────────────────────────────────────────────────────────────────────────────
# Inspect CSV markers + PDFs from Steps 6-7 then fill in below
# =============================================================================
cat("\n── Step 8: Annotation ──\n")
cat("Primate clusters:", levels(fib_primate$seurat_clusters), "\n")
cat("Rodent clusters: ", levels(fib_rodent$seurat_clusters),  "\n")

# ── Check Rat cell counts ─────────────────────────────────────────────────────
cat("\nRat subtype counts (from rodent object):\n")
print(table(fib_rodent$seurat_clusters[fib_rodent$species == "Rat"], useNA = "ifany"))

# ─── Fill in after marker inspection ─────────────────────────────────────────
primate_annotations <- c(
  "0"  = "Contaminant",
  "1"  = "Resident fibroblast",
  "2"  = "Adventitial fibroblast",
  "3"  = "Resident fibroblast",
  "4"  = "Contaminant",
  "5"  = "Contaminant",
  "6"  = "Contaminant",
  "7"  = "Resident fibroblast",
  "8"  = "Myofibroblast",
  "9"  = "Contaminant",
  "10" = "Contaminant",
  "11" = "Contaminant",
  "12" = "Contaminant",
  "13" = "Inflammatory CAF",
  "14" = "ECM-remodeling fibroblast",
  "15" = "Neural crest fibroblast",
  "16" = "Adventitial fibroblast",
  "17" = "ECM-remodeling fibroblast",
  "18" = "Contaminant"
)

rodent_annotations <- c(
  "0"  = "Inflammatory CAF",
  "1"  = "Contaminant",
  "2"  = "Contaminant",
  "3"  = "Resident fibroblast",
  "4"  = "ECM-remodeling fibroblast",
  "5"  = "Pericyte",
  "6"  = "Inflammatory CAF",
  "7"  = "Contaminant",
  "8"  = "Adventitial fibroblast",
  "9"  = "Contaminant",
  "10" = "Resident fibroblast",
  "11" = "ECM-remodeling fibroblast",
  "12" = "Contaminant",
  "13" = "Contaminant",
  "14" = "Adventitial fibroblast"
)

# Apply annotations
if (length(primate_annotations) > 0) {
  fib_primate$fib_subtype <- unname(
    primate_annotations[as.character(fib_primate$seurat_clusters)]
  )
} else {
  fib_primate$fib_subtype <- paste0("Cluster_", fib_primate$seurat_clusters)
  cat("  Primate: placeholder labels set — fill annotations above.\n")
}

if (length(rodent_annotations) > 0) {
  fib_rodent$fib_subtype <- unname(
    rodent_annotations[as.character(fib_rodent$seurat_clusters)]
  )
} else {
  fib_rodent$fib_subtype <- paste0("Cluster_", fib_rodent$seurat_clusters)
  cat("  Rodent: placeholder labels set — fill annotations above.\n")
}

# ── Broad category grouping ───────────────────────────────────────────────────
category_map <- c(
  "Pan-fibroblast"            = "Fibroblast",
  "Resident fibroblast"       = "Fibroblast",
  "Adventitial fibroblast"    = "Fibroblast",
  "ECM-remodeling fibroblast" = "Fibroblast",
  "Neural crest fibroblast"   = "Fibroblast",
  "Lipofibroblast"            = "Fibroblast",
  "Myofibroblast"             = "Activated stromal",
  "Inflammatory CAF"          = "Activated stromal",
  "Antigen-presenting CAF"    = "Activated stromal",
  "Pericyte"                  = "Mural",
  "Smooth muscle"             = "Mural",
  "Proliferating"             = "Cycling",
  "Contaminant"               = "Contaminant"
)

assign_category <- function(subtypes, cmap) {
  cat <- cmap[subtypes]
  cat[is.na(cat)] <- "Unannotated"
  unname(cat)
}

fib_primate$fib_category <- assign_category(fib_primate$fib_subtype, category_map)
fib_rodent$fib_category  <- assign_category(fib_rodent$fib_subtype,  category_map)

cat("\nPrimate fib_category distribution:\n")
print(table(fib_primate$fib_category, useNA = "ifany"))
cat("\nRodent fib_category distribution:\n")
print(table(fib_rodent$fib_category, useNA = "ifany"))

# =============================================================================
# STEP 9 — Annotated UMAPs (pre-integration)
# =============================================================================
cat("\n── Step 9: Annotated UMAPs ──\n")

plot_annotated <- function(seu, group_name, palette, spc_colors,
                           subtype_col = "fib_subtype") {
  p1 <- DimPlot(seu, group.by = subtype_col, cols = palette,
                label = TRUE, label.size = 3, repel = TRUE) +
    ggtitle(paste(group_name, "— fibroblast subtypes")) + theme_pub() +
    theme(legend.text = element_text(size = 8),
          legend.key.size = unit(0.4, "cm"))
  p2 <- DimPlot(seu, group.by = "species", cols = spc_colors) +
    ggtitle("Species") + theme_pub()
  p3 <- DimPlot(seu, group.by = "fib_category",
                label = TRUE, label.size = 3.5, repel = TRUE) +
    ggtitle("Category") + theme_pub()
  p1 | p2 | p3
}

cairo_pdf(file.path(outdir, "04_umap_primate_annotated.pdf"),
          width = 20, height = 6)
print(plot_annotated(fib_primate, "Primate (Human + Marmoset)",
                     subtype_palette, species_colors))
dev.off()

cairo_pdf(file.path(outdir, "04_umap_rodent_annotated.pdf"),
          width = 20, height = 6)
print(plot_annotated(fib_rodent, "Rodent (Mouse + Rat)",
                     subtype_palette, species_colors))
dev.off()

cat("  Annotated UMAPs saved.\n")

# =============================================================================
# STEP 10 — Merge primate + rodent → integrate all 4 species
# =============================================================================
cat("\n── Step 10: Merging primate + rodent for full integration ──\n")

fib_integrated <- merge(fib_primate, y = fib_rodent, merge.data = TRUE)
cat("Total cells after merge:", ncol(fib_integrated), "\n")

fib_integrated[["RNA"]] <- JoinLayers(fib_integrated[["RNA"]])
fib_integrated[["RNA"]] <- split(fib_integrated[["RNA"]],
                                 f = fib_integrated$species)

fib_integrated <- NormalizeData(fib_integrated, verbose = FALSE)
fib_integrated <- FindVariableFeatures(fib_integrated, nfeatures = 3000,
                                       verbose = FALSE)
fib_integrated <- ScaleData(fib_integrated, verbose = FALSE)
fib_integrated <- RunPCA(fib_integrated, npcs = N_PCS, verbose = FALSE)

cat("Running Harmony on merged object...\n")
fib_integrated <- IntegrateLayers(
  object         = fib_integrated,
  method         = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction  = "harmony",
  group.by.vars  = "species",
  verbose        = TRUE
)

fib_integrated[["RNA"]] <- JoinLayers(fib_integrated[["RNA"]])
fib_integrated <- FindNeighbors(fib_integrated, reduction = "harmony",
                                dims = 1:N_HARMONY_DIMS, verbose = FALSE)
fib_integrated <- FindClusters(fib_integrated, resolution = 0.5,
                               verbose = FALSE)
fib_integrated <- RunUMAP(fib_integrated, reduction = "harmony",
                          dims = 1:N_HARMONY_DIMS, verbose = FALSE)

cat("Integrated clusters:", length(levels(fib_integrated$seurat_clusters)), "\n")

# =============================================================================
# STEP 11 — Propagate annotations to integrated object
# =============================================================================
cat("\n── Step 11: Propagating annotations ──\n")

all_meta <- rbind(
  data.frame(cell         = colnames(fib_primate),
             fib_subtype  = fib_primate$fib_subtype,
             fib_category = fib_primate$fib_category,
             stringsAsFactors = FALSE),
  data.frame(cell         = colnames(fib_rodent),
             fib_subtype  = fib_rodent$fib_subtype,
             fib_category = fib_rodent$fib_category,
             stringsAsFactors = FALSE)
)

idx <- match(colnames(fib_integrated), all_meta$cell)
cat("Cells matched:", sum(!is.na(idx)), "/", ncol(fib_integrated), "\n")

fib_integrated$fib_subtype  <- all_meta$fib_subtype[idx]
fib_integrated$fib_category <- all_meta$fib_category[idx]

cat("Subtype distribution:\n")
print(table(fib_integrated$fib_subtype, useNA = "ifany"))

# =============================================================================
# STEP 12 — Annotated DotPlot per species
# =============================================================================
cat("\n── Step 12: Annotated DotPlots per species ──\n")

# Drop SCT/integrated from fib_integrated before dotplot
DefaultAssay(fib_integrated) <- "RNA"
if ("SCT" %in% Assays(fib_integrated))        fib_integrated[["SCT"]]        <- NULL
if ("integrated" %in% Assays(fib_integrated)) fib_integrated[["integrated"]] <- NULL
fib_integrated[["RNA"]] <- JoinLayers(fib_integrated[["RNA"]])

gene_annotation <- stack(subtype_markers) %>%
  rename(gene = values, subtype = ind) %>%
  mutate(subtype = as.character(subtype)) %>%
  group_by(gene) %>%
  slice(1) %>%
  ungroup() %>%
  filter(gene %in% genes_ordered)

gene_annotation$gene    <- factor(gene_annotation$gene,    levels = genes_ordered)
gene_annotation$subtype <- factor(gene_annotation$subtype, levels = names(subtype_markers))

species_list <- c("Human", "Marmoset", "Mouse", "Rat")

subtype_order_present <- intersect(subtype_order,
                                   unique(fib_integrated$fib_subtype[
                                     !is.na(fib_integrated$fib_subtype)]))

all_markers        <- unique(unlist(subtype_markers))
gene_presence      <- sapply(species_list, function(sp) {
  cells_sp <- colnames(fib_integrated)[fib_integrated$species == sp]
  rn       <- rownames(fib_integrated[, cells_sp])
  all_markers %in% rn
})
rownames(gene_presence) <- all_markers
genes_keep         <- all_markers[rowSums(gene_presence) >= 2]
genes_ordered_filt <- intersect(genes_ordered, genes_keep)

cat("\nGenes available for DotPlot:", length(genes_ordered_filt), "/",
    length(genes_ordered), "\n")

make_annotated_dotplot <- function(seu, species_name, genes,
                                   subtype_ord, ann_palette,
                                   gene_ann, subtype_col = "fib_subtype",
                                   base_size = 9) {
  cells_sp    <- colnames(seu)[seu$species == species_name]
  sub_obj     <- subset(seu, cells = cells_sp)
  genes_sp    <- intersect(genes, rownames(sub_obj))
  
  subtypes_sp <- intersect(subtype_ord,
                           unique(sub_obj[[subtype_col]][
                             !is.na(sub_obj[[subtype_col]])]))
  
  sub_obj$fib_subtype <- factor(sub_obj$fib_subtype, levels = subtypes_sp)
  Idents(sub_obj)     <- subtype_col
  
  cat(sprintf("  [%s] %d genes | %d subtypes\n",
              species_name, length(genes_sp), length(subtypes_sp)))
  
  dot_data <- suppressWarnings(
    DotPlot(sub_obj, features = genes_sp, scale = TRUE,
            col.min = -2, col.max = 2)$data
  ) %>%
    rename(gene = features.plot, subtype_y = id) %>%
    mutate(gene      = factor(gene, levels = genes_sp),
           subtype_y = factor(subtype_y, levels = subtypes_sp))
  
  p_main <- ggplot(dot_data,
                   aes(x = gene, y = subtype_y,
                       size = pct.exp, color = avg.exp.scaled)) +
    geom_point() +
    scale_color_gradient2(low = "#313695", mid = "white", high = "#a50026",
                          midpoint = 0, name = "Avg\nexp",
                          limits = c(-2, 2), oob = scales::squish) +
    scale_size_continuous(name = "% exp", range = c(0.2, 5),
                          breaks = c(10, 25, 50, 75, 100)) +
    labs(x = NULL, y = NULL, title = species_name) +
    theme_pub(base_size = base_size) +
    theme(axis.text.x      = element_text(angle = 45, hjust = 1,
                                          size = 6.5, face = "italic"),
          axis.text.y      = element_text(size = 8),
          panel.grid.major = element_line(color = "grey92", linewidth = 0.25),
          legend.key.size  = unit(0.3, "cm"),
          legend.text      = element_text(size = 6.5),
          legend.title     = element_text(size = 7),
          plot.title       = element_text(face = "bold", hjust = 0.5, size = 10),
          plot.margin      = margin(t = 2, r = 5, b = 2, l = 5))
  
  ann_data <- gene_ann %>%
    filter(gene %in% genes_sp) %>%
    filter(subtype %in% subtypes_sp) %>%
    mutate(gene = factor(gene, levels = genes_sp), y = 1)
  
  p_ann <- ggplot(ann_data, aes(x = gene, y = y, fill = subtype)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_manual(values = ann_palette, drop = FALSE) +
    scale_x_discrete(position = "top") +
    scale_y_continuous(expand = c(0, 0)) +
    theme_void() +
    theme(legend.position = "none",
          plot.margin = margin(t = 2, r = 5, b = 0, l = 5))
  
  label_data <- ann_data %>%
    filter(subtype %in% subtypes_sp) %>%
    mutate(gene_num = as.integer(gene)) %>%
    group_by(subtype) %>%
    summarise(xmid = mean(gene_num),
              xmin = min(gene_num) - 0.5,
              xmax = max(gene_num) + 0.5,
              .groups = "drop") %>%
    mutate(y = 1)
  
  p_label <- ggplot(label_data,
                    aes(x = xmid, y = y, label = subtype, color = subtype)) +
    geom_text(size = 2.2, fontface = "bold", hjust = 0.5, vjust = 0.5,
              show.legend = FALSE) +
    geom_segment(aes(x = xmin, xend = xmax, y = 0.45, yend = 0.45,
                     color = subtype), linewidth = 0.5, show.legend = FALSE) +
    scale_color_manual(values = ann_palette) +
    scale_x_continuous(limits = c(0.5, length(genes_sp) + 0.5), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1.5), expand = c(0, 0)) +
    theme_void() +
    theme(plot.margin = margin(t = 0, r = 5, b = 1, l = 5))
  
  p_ann / p_label / p_main + plot_layout(heights = c(0.06, 0.1, 1))
}

# Generate and save
fib_clean <- subset(fib_integrated,
                    subset = fib_subtype != "Contaminant")

plots_ann <- lapply(species_list, function(sp) {
  make_annotated_dotplot(
    seu          = fib_clean,
    species_name = sp,
    genes        = genes_ordered_filt,
    subtype_ord  = subtype_order_present,
    ann_palette  = subtype_palette,
    gene_ann     = gene_annotation
  )
})
names(plots_ann) <- species_list
plots_ann <- Filter(Negate(is.null), plots_ann)

for (sp in names(plots_ann)) {
  cairo_pdf(file.path(outdir, paste0("05_dotplot_annotated_", tolower(sp), ".pdf")),
            width = 18, height = 8)
  print(plots_ann[[sp]])
  dev.off()
  cat(sprintf("  Saved: 05_dotplot_annotated_%s.pdf\n", tolower(sp)))
}

cairo_pdf(file.path(outdir, "05_dotplot_annotated_all_species.pdf"),
          width = 20, height = 30)
print(wrap_plots(plots_ann, ncol = 1) +
        plot_annotation(
          title = "Fibroblast subtype marker expression per species",
          theme = theme(plot.title = element_text(face = "bold",
                                                  hjust = 0.5, size = 13))
        ))
dev.off()

# =============================================================================
# STEP 13 — Cell proportion barplot
# =============================================================================
cat("\n── Step 13: Cell proportion barplot ──\n")

prop_df <- fib_clean@meta.data %>%
  filter(!is.na(fib_subtype)) %>%
  group_by(species, fib_subtype) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(species) %>%
  mutate(total      = sum(n),
         proportion = n / total * 100) %>%
  ungroup() %>%
  mutate(fib_subtype = factor(fib_subtype, levels = subtype_order),
         species     = factor(species, levels = c("Human", "Marmoset", "Mouse", "Rat")))

totals_df <- prop_df %>% distinct(species, total)
label_df  <- prop_df %>%
  filter(proportion >= 5) %>%
  group_by(species) %>%
  arrange(fib_subtype) %>%
  mutate(label_y = cumsum(proportion) - proportion / 2) %>%
  ungroup()

p_bar <- ggplot(prop_df, aes(x = species, y = proportion, fill = fib_subtype)) +
  geom_bar(stat = "identity", width = 0.75, color = "white", linewidth = 0.2) +
  geom_text(data = label_df,
            aes(y = label_y, label = n),
            size = 2.8, color = "white", fontface = "bold") +
  geom_text(data = totals_df,
            aes(x = species, y = 103,
                label = paste0("n = ", comma(total)), fill = NULL),
            size = 3.2, fontface = "bold", color = "black") +
  scale_fill_manual(values = subtype_palette, name = "Fibroblast subtype") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.07)),
                     labels = function(x) paste0(x, "%")) +
  labs(x = NULL, y = "Cell proportion (%)",
       title = "Fibroblast subtype composition per species") +
  theme_pub(base_size = 11) +
  theme(legend.position    = "right",
        legend.text        = element_text(size = 9),
        legend.key.size    = unit(0.45, "cm"),
        axis.text.x        = element_text(size = 11, face = "bold"),
        panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3)) +
  guides(fill = guide_legend(ncol = 1))

cairo_pdf(file.path(outdir, "06_barplot_fib_proportion_by_species.pdf"),
          width = 10, height = 7)
print(p_bar)
dev.off()

write.csv(prop_df, file.path(outdir, "fib_subtype_counts_by_species.csv"),
          row.names = FALSE)
cat("  Barplot and count table saved.\n")

# =============================================================================
# STEP 14 — UMAP split by species
# =============================================================================
cat("\n── Step 14: UMAP split by species ──\n")

fib_clean$fib_subtype <- factor(fib_clean$fib_subtype, levels = subtype_order_present)
species_counts <- table(fib_clean$species)

panels <- lapply(c("Human", "Marmoset", "Mouse", "Rat"), function(sp) {
  n_cells  <- species_counts[sp]
  cells_sp <- colnames(fib_clean)[fib_clean$species == sp]
  DimPlot(fib_clean,
          group.by = "fib_subtype",
          cells    = cells_sp,
          cols     = subtype_palette,
          label    = FALSE, pt.size = 0.4) +
    ggtitle(paste0(sp, "\n(n = ", comma(n_cells), ")")) +
    theme_pub(base_size = 9) +
    theme(legend.position = "none",
          axis.text  = element_blank(),
          axis.ticks = element_blank(),
          axis.line  = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5, size = 10))
})

legend_plot <- DimPlot(fib_clean, group.by = "fib_subtype",
                       cols = subtype_palette) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3))) +
  theme(legend.text     = element_text(size = 9),
        legend.key.size = unit(0.45, "cm"))
shared_legend <- cowplot::get_legend(legend_plot)

cairo_pdf(file.path(outdir, "07_umap_split_species_fib_subtype.pdf"),
          width = 20, height = 5)
print((panels[[1]] | panels[[2]] | panels[[3]] | panels[[4]]) |
        wrap_elements(shared_legend))
dev.off()

# ── All-species combined UMAP ─────────────────────────────────────────────────
p_all <- DimPlot(fib_clean,
                 group.by  = "fib_subtype",
                 cols      = subtype_palette,
                 label     = TRUE, label.size = 3.5,
                 repel     = TRUE, pt.size = 0.3,
                 label.box = TRUE) +
  ggtitle("Fibroblast subtypes — all species") +
  theme_pub() +
  theme(legend.text     = element_text(size = 9),
        legend.key.size = unit(0.45, "cm"))

cairo_pdf(file.path(outdir, "08_umap_fib_subtypes_all.pdf"), width = 10, height = 8)
print(p_all)
dev.off()

cat("  UMAPs saved.\n")

# =============================================================================
# STEP 15 — Save all objects
# =============================================================================
cat("\n── Step 15: Saving objects ──\n")
saveRDS(fib_primate,    file.path(outdir, "fib_primate_annotated.rds"))
saveRDS(fib_rodent,     file.path(outdir, "fib_rodent_annotated.rds"))
saveRDS(fib_integrated, file.path(outdir, "fib_integrated_all4.rds"))
saveRDS(fib_clean,      file.path(outdir, "fib_clean_annotated.rds"))
cat("  All objects saved to:", outdir, "\n")

cat("\n── Session info ──\n")
sessionInfo()